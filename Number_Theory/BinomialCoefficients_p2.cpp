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
Time complexity O(p)
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

int count_ps_in_fact(ll n, ll p) {
  ll pp = p;
  int cnt = 0;
  while (pp <= n) {
    cnt += n / pp;
    pp = pp * p;
  }
  return cnt;
}

i128 p0pow[4 * 11234567];
void precompute(int p) {
  i128  = 1;
  ll pp = (ll)p * p;
  for (int i = 1; i < p; i++) {
    p0 = (i128)(p0)*i % pp;
  }
  p0pow[0] = 1;
  for (int i = 1; i < p; i++) {
    p0pow[i] = p0pow[i - 1] * p0 % pp;
    // cout<<i<<" :"<<(ll)p0pow[i]<<endl;
  }
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

ll factorial_stripped_p_mod_pp(ll n, ll p, string gap = "") {
  if (dp[n] > 0)
    return dp[n];
  tr_begin(dbg(n));
  if (n == 0)
    return 1;
  ll pp = (ll)p * p;
  // cout << gap + "factorial_stripped_p_mod_pp from " << n << "! p=" <<
  // pp
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
  ans =
      (i128)(ans)*factorial_stripped_p_mod_pp(n / p, p, gap + " ") % pp;
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
      ll first_residual_of_residual_no = lasti + p + 1;
      for (ll i = first_residual_of_residual_no; i <= n; i++) {
        // cout << "" << i << ",";
        ans = (i128)(ans)*i % pp;
      }
      // cout << endl;
      // if (ans==0)std::cout<<"Yo2 "<<ansb4<<"
      // "<<first_residual_of_residual_no<<"
      // "<<n<<std::endl;
    }
  }
  // cout << gap + "factorial_stripped_p_mod_pp from " << n << "! p=" <<
  // pp
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
  ll r_f_inv = inverse(factorial_stripped_p_mod_pp(r, p), p * p);
  // assert(1 <= r_f_inv);
  ll nr_f_inv = inverse(factorial_stripped_p_mod_pp(n - r, p), p * p);
  // assert(1 <= nr_f_inv);
  ll n_f = factorial_stripped_p_mod_pp(n, p);
  // assert(1 <= n_f);
  // cout << "inverse(factorial_stripped_p_mod_pp(" << r << ", " << p<<
  // "),
  // p*p) =" << r_f_inv << endl;
  // cout << "inverse(factorial_stripped_p_mod_pp(" << n - r << ", " <<
  // p<<
  // "), p*p) =" << nr_f_inv << endl;
  // cout << "factorial_stripped_p_mod_pp(" << n << ", " << p << ") ="
  // <<
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
  auto start = clock();

  // b = bonom_p2(240, 130, 7);
  // cout << "bonom_p2(240, 130) mod " << 7 * 7 << " = " << (ll)b << endl;
  // cout << "time from start= " << (double)(clock() - start) / CLOCKS_PER_SEC
  // << "s)" << endl;

  // p = 10007;
  // b = bonom_p2((ll)(1e17), (ll)(1e17) / 2, p);
  // cout << "bonom_p2(1e17, 5*1e16)%" << (ll)(p) * (p) << " = " << (ll)b
  // <<endl;
  // cout << "time from start= " << (double)(clock() - start) / CLOCKS_PER_SEC
  // << "s)" << endl;

  // p = 10007;
  // b = bonom_p2((ll)(99930048), (ll)(99930048) / 2, p);
  // cout << "bonom_p2(99930048, 99930048/2)%" << (ll)(p) * (p) << "=" <<(ll)b<<
  // endl;
  // cout << "time from start= " << (double)(clock() - start) / CLOCKS_PER_SEC
  // << "s)" << endl;

  // p = 100019;
  // b = bonom_p2((ll)(1e15), (ll)(1e15) / 2, p);
  // cout << "bonom_p2(1e15, 5*1e14)%"  << (ll)(p) * (p) << "=" << (ll)b <<endl;
  // cout << "time from start= " << (double)(clock() - start) / CLOCKS_PER_SEC
  // << "s)" << endl;

  p = 1000003;
  b = bonom_p2((ll)(1e14), (ll)(1e14) / 2, p);
  cout << "bonom_p2(1e14, 5*1e13)%" << (ll)(p) * (p) << "=" << (ll)b << endl;
  cout << "time from start= " << (double)(clock() - start) / CLOCKS_PER_SEC
       << "s)" << endl;

  p = 10000019;
  b = bonom_p2((ll)(1e17), (ll)(1e17) / 2, p);
  cout << "bonom_p2(1e17, 5*1e16)%" << (ll)(p) * (p) << "=" << (ll)b << endl;
  cout << "time from start= " << (double)(clock() - start) / CLOCKS_PER_SEC
       << "s)" << endl;
  return 0;
}
