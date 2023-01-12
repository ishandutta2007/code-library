#include <bits/stdc++.h>
using namespace std;

using ll = long long;
#define db long double
#define ii pair<int, int>
#define vi vector<ll>
#define fi first
#define se second
#define sz(a) (int)(a).size()
#define all(a) (a).begin(), (a).end()
#define pb push_back
#define mp make_pair
#define FN(i, n) for (int i = 0; i < (int)(n); ++i)
#define FEN(i, n) for (int i = 1; i <= (int)(n); ++i)
#define rep(i, a, b) for (int i = a; i < b; i++)
#define repv(i, a, b) for (int i = b - 1; i >= a; i--)
#define SET(A, val) memset(A, val, sizeof(A))

// const int N = 1e6 + 9;

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
vi rev;
ll FFTMOD;
void revbits(int newlim) {
  static int lim = -1;
  int t, j;
  if (newlim == lim)
    return;
  lim = newlim;
  rev.resize(lim + 1);
  int k = 0;
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

void fft(vector<cd> &poly, int inv = false) {
  int len, l;
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
        for (int j = 0; j < l; ++j) {
          if (j & 1)
            w[l + j] = w[(l + j) >> 1] * ww;
          else
            w[l + j] = w[(l + j) >> 1];
        }
      } else
        w[l] = cd(1.0, 0.0);
    }

    for (int i = 0; i < sz(poly); i += len)
      FN(j, l) {
        u = poly[i + j], v = poly[i + j + l] * w[l + j];
        poly[i + j] = u + v, poly[i + j + l] = u - v;
      }
  }
  if (inv)
    for (auto &x : poly)
      x = x / sz(poly);
}
vi multiply(vi &a, vi &b) {
  int bits = 1, sz1 = sz(a) + sz(b), reqsz;
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
  vi ans(sz1 - 1);
  FN(i, sz(ans)) ans[i] = (ll)(poly[i].real + 0.5) % FFTMOD;

  /*Uncomment for multiplication of two numbers
  int carry = 0;
  for (int i=0; i<(int)(ans.size()); ++i)
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
vi coeffs;
vi polynomial_chain_multiplication(int l, int r) {
  if (l == r) {
    vi tmp{-coeffs[l], 1};
    return tmp;
  }
  int mid = (l + r) >> 1;
  vi left = polynomial_chain_multiplication(l, mid);
  vi right = polynomial_chain_multiplication(mid + 1, r);
  vi ans = fft::multiply(left, right);
  return ans;
}
} // poly_chain



// void factorise(ll n) {}

int count_ps_in_fact(ll n, ll p) {
  ll pp = p;
  int cnt = 0;
  while (pp <= n) {
    cnt += n / pp;
    pp = pp * p;
  }
  return cnt;
}

// ll evaluate_P(Polynomial polyx, ll x, ll pp) {}

vi polyx;
void precompute(int p){
      fft::FFTMOD=(ll)(p)*p;
    // fft::FFTMOD=1e9+7;//(ll)(p)*p;
    poly_chain::coeffs.resize(0);
    for (int i = 0; i <= p - 1; i++)
      poly_chain::coeffs.pb(-i);
    polyx = poly_chain::polynomial_chain_multiplication(
        0, sz(poly_chain::coeffs) - 1);
    // cout << gap + "In " << n << "! " << pp << " Expanded: ";
    // for (int i = 0; i < polyx.size(); i++)
    //   cout << polyx[i] << " ";
    // cout << endl;
}

ll inverted_factorial_after_stripping(ll n, ll p, string gap = "") {
  if (n == 0)
    return 1;
  ll pp = (ll)pow(p, 2);
  cout << gap + "inverted_factorial_after_stripping from " << n << "! p=" << pp
       << endl;
  ll number_of_complete_groups = n / pp;
  ll residual_nos = n % pp;
  cout << gap + "In " << n << "! " << pp
       << " number_of_complete_groups =" << number_of_complete_groups
       << " residual_nos=" << residual_nos << endl;

  ll ans = (number_of_complete_groups % 2 == 0) ? 1 : -1;

  // group leaders ie k*PP as well as subleaders ie k*P are non-invertible so
  // handled resusively
  ans = ans * inverted_factorial_after_stripping(n / p, p, gap + " ") % pp;

  if (residual_nos <= p) {
    cout << gap + "In " << n << "! " << pp << " small residual_nos case from "
         << n - n % p + 1 << " to " << n << endl;
    // small residual_nos case
    for (ll i = n - n % p + 1; i <= n; i++) {
      ans = ans * i % pp;
    }
  } else {
    cout << gap + "In " << n << "! " << pp
         << " Large residual_nos case. Needs to be handled via polynomial"
         << endl;
    // Large residual_nos case.Needs to be handled via polynomial

    // P(x)=(x+1)(x+2)(x+3).....(x+p-1)
    cout << gap + "In " << n << "! " << pp
         << " We need to expand:(x+1)(x+2)...(x+" << p - 1 << ")" << endl;
    // note (x+p) should not be there as it is also not invertibe and should be
    // hadled recusively(already covered above)
    // Polynomial polyx
    // fft::FFTMOD=(ll)(p)*p;
    // // fft::FFTMOD=1e9+7;//(ll)(p)*p;
    // poly_chain::coeffs.resize(0);
    // for (int i = 0; i <= p - 1; i++)
    //   poly_chain::coeffs.pb(-i);
    // vi polyx = poly_chain::polynomial_chain_multiplication(
    //     0, sz(poly_chain::coeffs) - 1);
    cout << gap + "In " << n << "! " << pp << " Expanded: ";
    for (int i = 0; i < polyx.size(); i++)
      cout << polyx[i] << " ";
    cout << endl;

    // Evaluate polynomial at P(number_of_complete_groups*pp),
    // P(number_of_complete_groups*pp+p), P(number_of_complete_groups*pp+2*p)
    // upto P(number_of_complete_groups*k*pp)
    // such that number_of_complete_groups*k*pp+p-1<=n => k<=
    // (n+1-p)/(number_of_complete_groups*pp)
    cout << gap + "In " << n << "! " << pp << " We need to evaluate it at:";
    for (ll i = number_of_complete_groups * pp; i + p - 1 <= n; i += p) {
      cout << "P(" << i << "),";
      // ans = ans * evaluate_P(polyx, i, pp) % pp;
    }
    cout << endl;

    ll number_of_subgroups = residual_nos / p;
    cout << gap + "In " << n << "! " << pp << " small residual_nos stil left :";
    ll first_residual_no =
        number_of_complete_groups * pp + number_of_subgroups * p + 1;
    for (ll i = first_residual_no; i <= n; i++) {
      cout << "" << i << ",";
      ans = ans * i % pp;
    }
    cout << endl;
  }
  return ans;
}

ll bonom_pk(ll n, ll r, ll p, int k) {
  precompute(p);
  int cnt_n = count_ps_in_fact(n, p);
  cout << "cnt_n(" << n << "):" << cnt_n << endl;
  int cnt_r = count_ps_in_fact(r, p);
  cout << "cnt_r(" << r << "):" << cnt_r << endl;
  int cnt_nr = count_ps_in_fact(n - r, p);
  cout << "cnt_nr(" << n - r << "):" << cnt_nr << endl;
  int net_p = cnt_n - cnt_r - cnt_nr;
  if (net_p >= k)
    return 0;
  return inverted_factorial_after_stripping(r, p);
  // return 0;
}

// ll binom(ll n,ll r,ll m){
//   ans=1
//   map<ll, int> map_format =
//   FactorHelper::flat_to_map_format(PollardRho::factorize(m));
//   for (map<ll, int>::iterator it = map_format.begin(); it !=
//   map_format.end(); it++)
//       // cout << it->first << " : " << it->second << endl;
//       ans=ans*bonom_p^k(lln, llr, ll p,int k)
//   cout<<ans<<endl;
// }

int main() {
  // binom(10,10,1000);
  cout << bonom_pk(240, 130, 7, 2) << endl;
  return 0;
}
