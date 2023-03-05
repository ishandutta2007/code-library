#include <bits/stdc++.h>
/***
 * https://github.com/sgtlaugh/algovault/blob/master/code_library/fast_sieve.cpp
 *
 * Optimized Sieve of Eratosthenes
 *
 * Sieves all numbers from 1 to MAX and stores all the primes in primes[] array
 * The small primes are generated first using a simpler sieve
 * The numbers are chunked into blocks of fixed sizes
 * Each block is processed separately to make it cache-friendly
 * The is_composite[] array is a compressed bit-vector denoting the numbers
 *crossed out for each block
 * The sieve uses a wheel of size 15015 (3*5*7*11*13) to process each block
 *efficiently
 *
 * The algorithm can generate all the prime numbers from 1 to 2^31 in a little
 *under 1 seconds in a 4.00GHz core-i7 PC when compiled with -O2
 * Runtime in CodeForces - 1500 ms with GNU G++ 17
 *
 ***/

// takes 0.8s for n = 1e9

#define MAX 1000000010

using namespace std;

const uint32_t block_size = 1 << 19;

uint32_t s, prime_cnt, sq[65536], sp[65536], primes[2 * 26500010];
uint64_t wheel[15015], is_composite[8192], mask[12][62][8192];

inline void setbit(uint64_t *ar, uint32_t bit) {
  ar[bit >> 6] |= (1ULL << (bit & 63));
}

inline uint32_t get_idx(uint32_t i, uint32_t j) {
  if (sq[j] > i)
    return (sq[j] - i) >> 1;
  uint32_t x = sp[j] - i % sp[j];
  if (!(x & 1))
    x += sp[j];
  return x >> 1;
}

void small_sieve() {
  for (uint32_t i = 2; i * i < 65536; i++) {
    for (uint32_t j = i * i; j < 65536 && !sp[i]; j += i) {
      sp[j] = 1;
    }
  }
  for (uint32_t i = 2; i < 65536; i++) {
    if (!sp[i])
      sp[s] = i, sq[s++] = i * i;
  }
}

void process_block(uint32_t i) {
  uint32_t j, k, l, d, m, x, lim = i + block_size, idx = i % 15015, chunk = 0;

  idx = (idx + ((idx * 105) & 127) * 15015) >> 7;
  for (j = 0; (j << 7) < block_size; j += chunk, idx = 0) {
    chunk = min(15015 - idx, (block_size >> 7) - j);
    memcpy(is_composite + j, wheel + idx, sizeof(uint64_t) * chunk);
  }
  if (!i)
    is_composite[0] = (is_composite[0] | 1) & ~110ULL;

  l = block_size >> 1, m = block_size >> 7;
  for (j = 6; j < 18 && i; j++) {
    for (x = get_idx(i, j), k = 0, d = j - 6; k < m; k++) {
      is_composite[k] |= mask[d][x][k];
    }
  }

  for (j = (i == 0) ? 6 : 18; j < s && sq[j] < lim; j++) {
    for (x = get_idx(i, j); x < l; x += sp[j]) {
      setbit(is_composite, x);
    }
  }
}

void populate_primes(uint32_t i, uint32_t n) {
  for (uint32_t j = 0; (j << 7) < block_size; j++) {
    uint64_t x = ~is_composite[j];
    while (x) {
      uint32_t p = i + (j << 7) + (__builtin_ctzll(x) << 1) + 1;
      if (p <= n)
        primes[prime_cnt++] = p;
      x ^= (-x & x);
    }
  }
}

void fast_sieve(uint32_t n) {
  small_sieve();

  for (uint32_t i = 1; i <= 5; i++) {
    for (uint32_t j = i + (i > 3); j < 960960; j += sp[i]) {
      setbit(wheel, j);
    }
  }

  for (uint32_t i = 6; i <= 17; i++) {
    for (uint32_t j = 0; j < sp[i]; j++) {
      for (uint32_t k = j; k < (block_size >> 1); k += sp[i]) {
        setbit(mask[i - 6][j], k);
      }
    }
  }

  if (n >= 2)
    primes[prime_cnt++] = 2;
  for (uint32_t i = 0; i <= n; i += block_size) {
    process_block(i);
    populate_primes(i, n);
  }
}

namespace fio {
const int BUF_SIZE = 8192;

int buf_len = 0, inptr = 0, outptr = 0;
char inbuf[BUF_SIZE], outbuf[BUF_SIZE], tmpbuf[128];

inline char read_char() {
  if (inptr >= buf_len) {
    inptr = 0, buf_len = fread(inbuf, 1, BUF_SIZE, stdin);
    if (buf_len == 0)
      return EOF;
  }

  return inbuf[inptr++];
}

template <typename T,
          typename = typename enable_if<is_integral<T>::value, T>::type>
bool read_one(T &x) {
  int c = ' ', neg = 0;
  while (c != '-' && !isdigit(c) && c != EOF)
    c = read_char();
  if (c == '-')
    neg = 1, c = read_char();
  if (c == EOF)
    return false;

  for (x = 0; isdigit(c); c = read_char()) {
    x = x * 10 + c - '0';
  }
  if (neg)
    x = -x;

  return true;
}

bool read_one(string &s) {
  int c = ' ';
  while (isspace(c) && c != EOF)
    c = read_char();
  if (c == EOF)
    return false;

  for (s.clear(); !isspace(c); c = read_char()) {
    s.push_back(c);
  }
  return true;
}

template <typename T> bool read_one(vector<T> &v, bool read_length = true) {
  if (read_length) {
    int n;
    if (!read_one(n))
      return false;
    v.resize(n);
  }

  for (auto &&x : v) {
    if (!read_one(x))
      return false;
  }
  return true;
}

bool read_line(string &s) {
  s.clear();
  int c = '\n';
  while (c == '\n' || c == '\r')
    c = read_char();

  while (c != '\n' && c != '\r' && c != EOF) {
    s.push_back(c);
    c = read_char();
  }
  return c != EOF;
}

int read() { return 0; }

template <typename T, typename... Args> int read(T &x, Args &...args) {
  if (!read_one(x))
    return 0;
  return read(args...) + 1;
}

/* End of read methods, write methods start below */

void flush() {
  fwrite(outbuf, 1, outptr, stdout);
  outptr = 0;
}

inline void write_char(const char &c) {
  if (outptr == BUF_SIZE)
    flush();
  outbuf[outptr++] = c;
}

void write_one(const char *s) {
  for (int j = 0; s[j]; j++)
    write_char(s[j]);
}

void write_one(const string &s) {
  for (auto &&c : s)
    write_char(c);
}

template <typename T,
          typename = typename enable_if<is_integral<T>::value, T>::type>
void write_one(T x) {
  if (x < 0)
    x = -x, write_char('-');

  int l = 0;
  while (x || !l) {
    tmpbuf[l++] = (x % 10) + '0';
    x /= 10;
  }
  while (l)
    write_char(tmpbuf[--l]);
}

template <typename T> void write_one(const vector<T> &v) {
  for (int i = 0; i < (int)v.size(); i++) {
    if (i)
      write_char(' ');
    write_one(v[i]);
  }
}

void write() {}

template <typename T, typename... Args>
void write(const T x, const Args... args) {
  write_one(x);
  write_char(sizeof...(args) && is_trivial<T>::value ? ' ' : '\n');
  write(args...);
}
} // namespace fio

int main() {
  int n = 999999999;
  fast_sieve(n);
  return 0;
}
