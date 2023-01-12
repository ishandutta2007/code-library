#include <bits/stdc++.h>
using namespace std;

const int N = 1e6 + 9;
using ll = long long;

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
