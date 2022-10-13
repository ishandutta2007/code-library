// nc0 + nc1 + nc2 + ... + nck
ll sum_bin_coef(ll n, ll k) { // O(log (k))
  ll i, bincoef = 1, sum = 1;
  for (i = 1; i <= k; ++i) {
    bincoef = bincoef * (n - i + 1) / i;
    sum += bincoef;
  }
  // cout << "sum_bin_coef @(" << k << "," << n << "):" << sum << endl;
  return sum;
}
