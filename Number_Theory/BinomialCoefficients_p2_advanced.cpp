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
Time complexity O(sqrt(p)log^2(p))
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

i128 fact0pow[4 * 11234567];
vector<i128> pans1_sqrtblocks;
vector<i128> pans2_sqrtblocks;

//O(P)
void precompute(int p) {
  auto precompute_start = clock();

  ll pp = (ll)p * p;
  i128 fact0 = 1;
  for (int i = 1; i < p; i++) {
    fact0 = (i128)(fact0)*i % pp;
  }
  fact0pow[0] = 1;
  for (int i = 1; i < p; i++) {
    fact0pow[i] = fact0pow[i - 1] * fact0 % pp;
    // cout<<i<<" :"<<(ll)fact0pow[i]<<endl;
  }
  int sqrtp=(int)sqrt(p);
  int block_size = sqrtp;

  ll sumofinverses = 0;
  i128 pans2 = 1;
  for(int k=1;k<=p/block_size;k++){
    for(int i=1;i<=block_size;i++){
      ll next_term=block_size * (k - 1) + i;
      pans2 = pans2 * next_term % pp;
      sumofinverses = (sumofinverses + inverse(next_term,  pp)) % pp;
    }
    i128 pans1 = (i128)p;
    pans1 = pans1 * pans2 % pp;
    pans1 = pans1 * (i128)sumofinverses % pp;
    // cout<< "of prime "<<p<<" precomptued "<<k<<"-th sqrt block prods="<< (ll)pans1<<" "<<(ll)pans2<<endl;
    pans1_sqrtblocks.push_back(pans1);
    pans2_sqrtblocks.push_back(pans2);
  }

  cout << "precompute time= " << (double)(clock() - precompute_start) / CLOCKS_PER_SEC  << "s)" << endl;
}

//T(N)=O(sqrtP)+T(N/P)
//T(N/P)=O(sqrtP)+T(N/PP)
//T(N/PP)=O(sqrtP)+T(N/PPP)
//O(sqrtP*log(N/P)/log(P))
map<ll, ll> dp;
ll factorial_stripped_p_mod_pp(ll n, ll p, string gap = "") {
  if (dp[n] > 0)
    return dp[n];
  // tr_begin(dbg(n));
  if (n == 0)
    return 1;
  ll pp = (ll)p * p;
  int sqrtp = (int)sqrt(p);
  ll u = n / pp; // number_of_complete_groups
  ll v = (n % pp) / p; // number_of_complete_subgroups
  ll w = (n % p) / sqrtp; // number_of_sqrtp_blocks
  ll r = ((n % pp) % p) % sqrtp; // number_of_final_residual numbers

  i128 ans = (u % 2 == 0) ? (i128)1 : (i128)(pp - 1);
  // cout << gap + "In " << n << "! no_of_groups=" << u << ", "<< "ans after group="<< (ll)ans << endl;
  ans = ans*fact0pow[v] % pp;
  // cout << gap + "In " << n << "! no_of_sungroups=" << v << ", "<< "ans after subgroup="<< (ll)ans << endl;
  if (w>=1) ans = ans*((i128)v * pans1_sqrtblocks[w-1] % pp + pans2_sqrtblocks[w-1]) % pp;
  // cout << gap + "In " << n << "! no_of_sqrtp_blocks=" << w << ", "<< v <<"*"<< (ll)pans1_sqrtblocks[w-1] <<" + "<< (ll)pans2_sqrtblocks[w-1] << endl;
  // cout << gap + "In " << n << "! " << "ans after sqrtp blocks="<< (ll)ans << endl;
  if (r>=1)
  for (ll i = n - r + 1; i <= n; i++) {
    ans = ans*(i128)i % pp;
  }
  // cout << gap + "In " << n << "! " << "ans after residual="<< (ll)ans << endl;
  if(n >= p)
  ans = ans*(i128)factorial_stripped_p_mod_pp(n / p, p, gap + " ") % pp;
  // cout << gap + "In " << n << "! " << "ans final after recusrive="<< (ll)ans << endl;
  dp[n] = (ll)ans;
  return (ll)ans;
  // return tr_end(ans, dbg(n));
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

  // p = 1000003;
  // b = bonom_p2((ll)(1e14), (ll)(1e14) / 2, p);
  // cout << "bonom_p2(1e14, 5*1e13)%" << (ll)(p) * (p) << "=" << (ll)b <<endl;
  // cout << "time from start= " << (double)(clock() - start) / CLOCKS_PER_SEC
  // << "s)" << endl;

  p = 10000019;
  b = bonom_p2((ll)(1e17), (ll)(1e17) / 2, p);
  cout << "bonom_p2(1e17, 5*1e16)%" << (ll)(p) * (p) << "=" << (ll)b << endl;
  cout << "time from start= " << (double)(clock() - start) / CLOCKS_PER_SEC
       << "s)" << endl;
  return 0;
}
//https://web.archive.org/web/20170202003812/http://www.dms.umontreal.ca/~andrew/PDF/BinCoeff.pdf
//https://blog.prabowodjonatan.id/posts/binomial-mod-pe/
//https://en.wikipedia.org/wiki/Wilf%E2%80%93Zeilberger_pair

