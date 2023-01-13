#define TRACE
#ifdef TRACE
#include <iostream>
#include <sstream>
class tracer_t {
  int trace_depth = 0;
  inline void prefix(int depth, int line, const std::string &fun) {
    std::cout << std::string(depth, ' ') << depth << ":LINE " << line << ": "
              << fun << "(";
  }
  inline void suffix() { std::cout << ") "; }
  template <class Arg, class... Args>
  inline void suffix(const Arg &arg, const Args &... args) {
    if (std::cout << stringify(arg), sizeof...(Args) > 0)
      std::cout << ",";
    suffix(args...);
  }
  inline void pause() { std::cout << std::endl, std::cin.get(); }

public:
  template <class Arg> inline std::string stringify(const Arg &arg) {
    std::stringstream ss;
    ss << arg;
    const auto s = ss.str();
    return std::is_same<Arg, std::string>::value ? "\"" + s + "\"" : s;
  }
  template <class... Args>
  inline void begin(int line, const std::string &fun, const Args &... args) {
    prefix(++trace_depth, line, fun), suffix(args...), std::cout << "called",
        pause();
  }
  template <class... Args>
  inline void trace(int line, const std::string &fun, const Args &... args) {
    prefix(trace_depth, line, fun), suffix(args...), pause();
  }
  template <class Ans, class... Args>
  inline Ans end(int line, const std::string &fun, const Ans &ans,
                 const Args &... args) {
    prefix(trace_depth--, line, fun), suffix(args...),
        std::cout << "return value = " << ans, pause();
    return ans;
  }
} tracer;
#define dbg(arg) string(#arg) + " = " + tracer.stringify(arg)
#define tr_begin(...) tracer.begin(__LINE__, __FUNCTION__, __VA_ARGS__)
#define tr(...) tracer.trace(__LINE__, __FUNCTION__, __VA_ARGS__)
#define tr_end(...) tracer.end(__LINE__, __FUNCTION__, __VA_ARGS__)
#else
#define tr_begin(...)
#define tr(...)
#define tr_end(ans, ...) ans
#define dbg(arg) arg
#endif // TRACE
/*
ncr modulo p^2 where p is 10^7
n and r are also 64 bit
*/
#include <bits/stdc++.h>
using namespace std;

using i128 = __int128_t;
using ll = long long;
using db = long double;
#define ii pair<int, int>
#define vll vector<ll>
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
vll rev;
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
vll multiply(vll &a, vll &b) {
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
  vll ans(sz1 - 1);
  FN(i, sz(ans)) ans[i] = (ll)((i128)(poly[i].real + 0.5) % FFTMOD);

  /*Uncomment for multiplication of two numbers
  int carry = 0;
  for (int i=0; i<(int)(ans.size()); ++i)
  {
    ans[i] += carry;
    carry = ans[i] / 10;
    ans[i] %= 10;
  }
  */
  bool hasNegative =
      std::any_of(ans.begin(), ans.end(), [](int x) { return x < 0; });
  if (hasNegative) {
    cout << "nega" << a.size() << " " << b.size() << endl;
    for (int i = 0; i < a.size(); i++)
      cout << a[i] << " ";
    cout << endl;
    for (int i = 0; i < b.size(); i++)
      cout << b[i] << " ";
    cout << endl;
    for (int i = 0; i < ans.size(); i++)
      cout << poly[i].real << " " << ans[i] << " ,";
    cout << endl;
  }

  return ans;
}
} // fft

namespace poly_chain {
vll coeffs;
vll polynomial_chain_multiplication(int l, int r) {
  if (l == r) {
    vll tmp{coeffs[l], 1};
    return tmp;
  }
  int mid = (l + r) >> 1;
  vll left = polynomial_chain_multiplication(l, mid);
  vll right = polynomial_chain_multiplication(mid + 1, r);
  vll ans = fft::multiply(left, right);
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

i128 p0;
// vll polyx;
i128 p0pow[4 * 11234567];
void precompute(int p) {
  p0 = 1;
  ll pp = (ll)p * p;
  for (int i = 1; i < p; i++) {
    p0 = (i128)(p0)*i % pp;
  }
  p0pow[0] = 1;
  for (int i = 1; i < p; i++) {
    p0pow[i] = p0pow[i - 1] * p0 % pp;
    // cout<<i<<" :"<<(ll)p0pow[i]<<endl;
  }
  // fft::FFTMOD = (ll)(p)*p;
  // // fft::FFTMOD=1e9+7;//(ll)(p)*p;
  // poly_chain::coeffs.resize(0);
  // for (int i = 0; i <= p - 1; i++)
  //   poly_chain::coeffs.pb(i);
  // polyx = poly_chain::polynomial_chain_multiplication(
  //     0, sz(poly_chain::coeffs) - 1);
  // for (int i = 0; i < polyx.size(); i++)
  //   cout << polyx[i] << " ";
  // cout << endl;
  // powsp[0]=1;
  // for(int i=1;i<p;i++)
  // powsp[i]=(powsp[i-1]*p)%pp;
}

ll bigpow(ll n, ll k, ll MOD) {
  ll x = 1, y = k;
  while (y) {
    if (y & 1) {
      x = x * n;
      if (x >= MOD)
        x = x % MOD;
    }
    n = n * n;
    if (n >= MOD)
      n = n % MOD;
    y >>= 1;
  }
  return x;
}

std::map<ll, ll> dp;

ll factorial_after_stripping_mod_p2(ll n, ll p, string gap = "") {
  if (dp[n] > 0)
    return dp[n];
  tr_begin(dbg(n));
  if (n == 0)
    return 1;
  ll pp = (ll)p * p;
  // cout << gap + "factorial_after_stripping_mod_p2 from " << n << "! p=" << pp
  // << " starts" << endl;
  ll number_of_complete_groups = n / pp;
  ll residual_nos = n % pp;
  // cout << gap + "In " << n << "! " << pp
  //      << " number_of_complete_groups =" << number_of_complete_groups
  //      << " residual_nos=" << residual_nos << endl;

  ll ans = (number_of_complete_groups % 2 == 0) ? 1 : pp - 1;

  // group leaders ie k*PP as well as subleaders ie k*P are non-invertible so
  // handled resusively
  // ll ansb4=ans;
  ans = (i128)(ans)*factorial_after_stripping_mod_p2(n / p, p, gap + " ") % pp;
  // if (ans==0)std::cout<<"Yo0"<<ansb4<<" "<<n<<std::endl;
  // cout << gap + "In " << n << "! " << pp << " factorial from recursive case
  // of "
  //      << n << " ie " << n / p << " is " << ans << endl;
  if (residual_nos <= p) {
    // cout << gap + "In " << n << "! " << pp << " small residual_nos case from
    // "
    //      << n - n % p + 1 << " to " << n << endl;
    // small residual_nos case
    for (ll i = n - n % p + 1; i <= n; i++) {
      ans = (i128)(ans)*i % pp;
    }
  } else {
    // cout << gap + "In " << n << "! " << pp
    //      << " Large residual_nos case. Needs to be handled via polynomial"
    //      << endl;
    // Large residual_nos case.Needs to be handled vlla polynomial

    // P(x)=(x+1)(x+2)(x+3).....(x+p-1)
    // cout << gap + "In " << n << "! " << pp
    //      << " We need to expand:(x+1)(x+2)...(x+" << p - 1 << ")" << endl;
    // note (x+p) should not be there as it is also not invertibe and should be
    // hadled recusively(already covered above)
    // Polynomial polyx
    // fft::FFTMOD=(ll)(p)*p;
    // // fft::FFTMOD=1e9+7;//(ll)(p)*p;
    // poly_chain::coeffs.resize(0);
    // for (int i = 0; i <= p - 1; i++)
    //   poly_chain::coeffs.pb(-i);
    // vll polyx = poly_chain::polynomial_chain_multiplication(
    //     0, sz(poly_chain::coeffs) - 1);
    // cout << gap + "In " << n << "! " << pp << " Expanded: ";
    // for (int i = 1; i < polyx.size(); i++)
    //   cout << polyx[i] << " ";
    // cout << endl;
    // ll p0 = polyx[1];
    // ll p0 = 1;
    // for (int i = 1; i < p; i++)
    //   p0 = (i128)(p0)*i % pp;

    // Evaluate polynomial at P(number_of_complete_groups*pp),
    // P(number_of_complete_groups*pp+p), P(number_of_complete_groups*pp+2*p)
    // upto P(number_of_complete_groups*k*pp)
    // such that number_of_complete_groups*k*pp+p-1<=n => k<=
    // (n+1-p)/(number_of_complete_groups*pp)
    // cout << gap + "In " << n << "! " << pp << " We need to evaluate it at:";

    ll lasti = -1;
    // ll actual = 0;
    // ll ansb4=ans;
    // for (ll i = number_of_complete_groups * pp; i + p - 1 <= n; i += p) {
    //   // cout << "P(" << i << "),";
    // lasti = i;
    //   actual++;
    //   ans = (i128)(ans)*p0 % pp;
    // }
    // cout << endl;
    ll p0times = ((n + 1 - p) - (number_of_complete_groups * pp)) / p + 1;
    // cout << "p0times=" << p0times << " actual=" << actual << " n=" << n <<
    // endl;
    // assert(p0times == actual);
    // cout<<ans<<" "<<ansb4<<"*"<<(ll)p0pow[p0times]<<"
    // would="<<(ll)((i128)ansb4*p0pow[p0times] % pp)<<endl;
    if (p0times >= 1) {
      ans = (i128)(ans)*p0pow[p0times] % pp;
      //   ans = (i128)(ans)*bigpow(p0, p0times, pp) % pp;
      lasti = number_of_complete_groups * pp + (p0times - 1) * p;
    }
    // if (ans==0)std::cout<<"Yo1"<<std::endl;

    ll number_of_subgroups = residual_nos / p;
    if (lasti >= 0) {
      // ll ansb4=ans;
      // cout << gap + "In " << n << "! " << pp << " small residual_nos stil
      // left :";
      ll first_residual_no = lasti + p + 1;
      for (ll i = first_residual_no; i <= n; i++) {
        // cout << "" << i << ",";
        ans = (i128)(ans)*i % pp;
      }
      // cout << endl;
      // if (ans==0)std::cout<<"Yo2 "<<ansb4<<" "<<first_residual_no<<"
      // "<<n<<std::endl;
    }
  }
  // cout << gap + "factorial_after_stripping_mod_p2 from " << n << "! p=" << pp
  //      << " ends with return value=" << ans << endl;

  // return ans;
  dp[n] = ans;
  return tr_end(ans, db(n));
}

ll extended_euclid(ll a, ll b, ll &x, ll &y) {
  if (b == 0) {
    x = 1;
    y = 0;
    return a;
  }
  ll x1, y1;
  ll d = extended_euclid(b, a % b, x1, y1);
  x = y1;
  y = x1 - y1 * (a / b);
  return d;
}

ll inverse(ll a, ll m) {
  ll x, y;
  ll g = extended_euclid(a, m, x, y);
  if (g != 1)
    return -1;
  return (x % m + m) % m;
}

i128 bonom_p2(ll n, ll r, ll p) {
  precompute(p);
  ll pp = (ll)p * p;
  int cnt_n = count_ps_in_fact(n, p);
  // cout << "cnt_n(" << n << "):" << cnt_n << endl;
  int cnt_r = count_ps_in_fact(r, p);
  // cout << "cnt_r(" << r << "):" << cnt_r << endl;
  int cnt_nr = count_ps_in_fact(n - r, p);
  // cout << "cnt_nr(" << n - r << "):" << cnt_nr << endl;
  int net_p = cnt_n - cnt_r - cnt_nr;
  // cout << "Excess p on numeratot=" << net_p << endl;
  if (net_p >= 2)
    return 0;
  ll r_f_inv = inverse(factorial_after_stripping_mod_p2(r, p), p * p);
  // assert(1 <= r_f_inv);
  ll nr_f_inv = inverse(factorial_after_stripping_mod_p2(n - r, p), p * p);
  // assert(1 <= nr_f_inv);
  ll n_f = factorial_after_stripping_mod_p2(n, p);
  // assert(1 <= n_f);
  // cout << "inverse(factorial_after_stripping_mod_p2(" << r << ", " << p<< "),
  // p*p) =" << r_f_inv << endl;
  // cout << "inverse(factorial_after_stripping_mod_p2(" << n - r << ", " << p<<
  // "), p*p) =" << nr_f_inv << endl;
  // cout << "factorial_after_stripping_mod_p2(" << n << ", " << p << ") =" <<
  // n_f<< endl;
  i128 ans = n_f;
  ans = ans * r_f_inv % pp;
  ans = ans * nr_f_inv % pp;
  // cout<<ans;
  return (net_p == 1) ? ans * p % pp : ans;
  // return 0;
}

int main() {
  ll p;
  i128 b;
  // b = bonom_p2(240, 130, 7);
  // cout << "bonom_p2(240, 130) mod " << 7 * 7 << " = " << (ll)b << endl;

  // p = 10007;
  // b = bonom_p2((ll)(1e17), (ll)(1e17) / 2, p);
  // cout << "bonom_p2(1e17, 5*1e16)%" << (ll)(p) * (p) << " = " << (ll)b
  // <<endl;

  // p = 10007;
  // b = bonom_p2((ll)(99930048), (ll)(99930048) / 2, p);
  // cout << "bonom_p2(99930048, 99930048/2)%" << (ll)(p) * (p) << "=" <<(ll)b<<
  // endl;

  // p = 100019;
  // b = bonom_p2((ll)(1e15), (ll)(1e15) / 2, p);
  // cout << "bonom_p2(1e15, 5*1e14)%"  << (ll)(p) * (p) << "=" << (ll)b <<endl;

  // p = 1000003;
  // b = bonom_p2((ll)(1e14), (ll)(1e14) / 2, p);
  // cout << "bonom_p2(1e14, 5*1e13)%" << (ll)(p) * (p) << "=" << (ll)b <<endl;

  p = 10000019;
  b = bonom_p2((ll)(1e17), (ll)(1e17) / 2, p);
  cout << "bonom_p2(1e17, 5*1e16)%" << (ll)(p) * (p) << "=" << (ll)b << endl;
  return 0;
}
