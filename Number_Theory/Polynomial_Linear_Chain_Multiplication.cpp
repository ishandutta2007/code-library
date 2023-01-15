/*
 Polynomial Chain Multiplication
 Takes input n
 (also kept provistion to input vector<int>coeffs -> coefficients of
 polynomials)
 Computes (1+x)(2+x)...(n+x)
 Returns vector<ll>ans -> coefficients of resultant polynomial
 O(N*log^2(N))
*/

#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace std;
using namespace __gnu_pbds;

using i128 = __int128_t;
using ll = long long;
#define db long double
#define vi128 vector<i128>
#define sz(a) (int)(a).size()
#define pb push_back
#define FN(i, n) for (int i = 0; i < (int)(n); ++i)
#define FEN(i, n) for (int i = 1; i <= (int)(n); ++i)

const db PI = acos((db)-1);
namespace fft {
template <class T> class comp {
public:
  T real, img;
  comp(T a = (T)0, T b = (T)0) : real(a), img(b) {}

  comp conj() { return comp(this->real, -(this->img)); }
  comp operator=(const comp &a) {
    this->real = a.real, this->img = a.img;
    return *this;
  }
  comp operator+(const comp &b) {
    return comp(this->real + b.real, this->img + b.img);
  }
  comp operator-(const comp &b) {
    return comp(this->real - b.real, this->img - b.img);
  }
  comp operator*(const T &num) {
    return comp(this->real * num, this->img * num);
  }
  comp operator/(const T &num) {
    return comp(this->real / num, this->img / num);
  }
  comp operator*(const comp &b) {
    return comp(this->real * b.real - this->img * b.img,
                this->img * b.real + this->real * b.img);
  }
  comp operator/(const comp &b) {
    comp temp(b.real, -b.img);
    comp n = (*this) * temp;
    return n / (b.x * b.x + b.y * b.y);
  }
};
#define cd comp<double>
vector<cd> w;
vi128 rev;
i128 FFTMOD;
void revbits(i128 newlim) {
  static i128 lim = -1;
  i128 t, j;
  if (newlim == lim)
    return;
  lim = newlim;
  rev.resize(lim + 1);
  i128 k = 0;
  while ((1 << k) < newlim)
    ++k;
  assert((1 << k) == newlim);
  FEN(i, lim) {
    j = rev[i - 1];
    t = k - 1;
    while (t >= 0 && ((j >> t) & 1))
      j ^= (1 << t), --t;
    if (t >= 0)
      j ^= (1 << t), --t;
    rev[i] = j;
  }
}

void fft(vector<cd> &poly, i128 inv = false) {
  i128 len, l;
  revbits(sz(poly));
  if (inv)
    for (auto &x : poly)
      x = x.conj();

  FN(i, sz(poly)) if (rev[i] > i) swap(poly[i], poly[rev[i]]);
  cd u, v;
  if (sz(w) < sz(poly))
    w.resize(sz(poly));
  for (len = 2, l = 1; len <= sz(poly); len += len, l += l) {
    if (w[l].real == 0 && w[l].img == 0) {
      double ang = PI / l;
      cd ww(cos(ang), sin(ang));
      if (l > 1) {
        for (i128 j = 0; j < l; ++j) {
          if (j & 1)
            w[l + j] = w[(l + j) >> 1] * ww;
          else
            w[l + j] = w[(l + j) >> 1];
        }
      } else
        w[l] = cd(1.0, 0.0);
    }

    for (i128 i = 0; i < sz(poly); i += len)
      FN(j, l) {
        u = poly[i + j], v = poly[i + j + l] * w[l + j];
        poly[i + j] = u + v, poly[i + j + l] = u - v;
      }
  }
  if (inv)
    for (auto &x : poly)
      x = x / sz(poly);
}
vi128 multiply(vi128 &a, vi128 &b) {
  i128 bits = 1, sz1 = sz(a) + sz(b), reqsz;
  while ((1 << bits) < sz1)
    ++bits;
  reqsz = (1 << bits);
  vector<cd> poly(reqsz);
  FN(i, sz(a)) poly[i].real = a[i];
  FN(i, sz(b)) poly[i].img = b[i];
  fft(poly);
  cd p, qtmp, q;
  poly[0] = poly[0].real * poly[0].img; // for i = 0
  FEN(i, reqsz >> 1) {
    p = poly[i] + poly[reqsz - i].conj(),
    qtmp = poly[reqsz - i] - poly[i].conj();
    q.real = qtmp.img, q.img = qtmp.real;
    poly[i] = (p * q) * 0.25;
    if (i)
      poly[reqsz - i] = poly[i].conj();
  }
  fft(poly, true);
  vi128 ans(sz1 - 1);
  FN(i, sz(ans)) ans[i] = (i128)(poly[i].real + 0.5) % FFTMOD;

  /*Uncomment for multiplication of two numbers
  i128 carry = 0;
  for (i128 i=0; i<(i128)(ans.size()); ++i)
  {
    ans[i] += carry;
    carry = ans[i] / 10;
    ans[i] %= 10;
  }
  */
  return ans;
}
} // fft

namespace poly_chain {
vi128 coeffs;
vi128 polynomial_chain_multiplication(int l, int r) {
  if (l == r) {
    vi128 tmp{coeffs[l], 1};
    return tmp;
  }
  i128 mid = (l + r) >> 1;
  vi128 left = polynomial_chain_multiplication(l, mid);
  vi128 right = polynomial_chain_multiplication(mid + 1, r);
  vi128 ans = fft::multiply(left, right);
  return ans;
}
} // poly_chain

int32_t main() {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  cout.tie(NULL);
  int n;
  cin >> n;
  for (int i = 0; i < n; i++)
    poly_chain::coeffs.pb(i + 1);
  fft::FFTMOD = 100000380000361L; //(1e7+19)*(1e7+19);
  vi128 ans = poly_chain::polynomial_chain_multiplication(
      0, sz(poly_chain::coeffs) - 1);

  for (int i = 0; i < ans.size(); i++) {
    cout << (ll)ans[i] << "x^" << i;
    if (i < ans.size() - 1)
      cout << " + ";
  }
  cout << endl;
  return 0;
}
