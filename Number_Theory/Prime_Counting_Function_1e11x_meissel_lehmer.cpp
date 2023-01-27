#include <bits/stdc++.h>
using namespace std;
#define int long long
typedef long long ll;
 
const int N = 1911570 + 7, M = 632454 + 7, K = 10, P = 2 * 3 * 5 * 7 * 11 * 13, Q = 6, R = 29;
const double eps = 1e-18;
ll sqrt_n, limit;
int pi_list[N], id1[M], id2[M], pre_phi[K + 7][P + 7], product[Q + 7];
ll prime[N], number[M], big_pi_list[M];
double inv[N];
bool p[N], vis[M];
 
inline ll sqrt(ll n){
	ll ans = sqrt((double)n);
	while (ans * ans <= n) ans++;
	return ans - 1;
}
 
inline ll max(ll a, ll b){
	return a > b ? a : b;
}
 
inline void init(ll n){
	int cnt = 0, id = 0;
	sqrt_n = sqrt(n);
	limit = max(sqrt_n * sqrt(log2(n)), R);
	p[0] = p[1] = true;
	pi_list[1] = 0;
	for (register ll i = 2; i <= limit; i++){
		if (!p[i]){
			cnt++;
			prime[cnt] = i;
			inv[cnt] = 1.0 / i + eps;
		}
		pi_list[i] = cnt;
		for (register int j = 1; j <= cnt && i * prime[j] <= limit; j++){
			p[i * prime[j]] = true;
			if (i % prime[j] == 0) break;
		}
	}
	for (register ll i = 1, j; i <= n; i = j + 1){
		ll tn = n / i;
		j = n / tn;
		id++;
		number[id] = tn;
		if (tn <= sqrt_n){
			id1[tn] = id;
		} else {
			id2[j] = id;
		}
	}
	for (register int i = 1; i <= P; i++){
		pre_phi[0][i] = i;
	}
	for (register int i = 1; i <= K; i++){
		for (register int j = 1; j <= P; j++){
			pre_phi[i][j] = pre_phi[i - 1][j] - pre_phi[i - 1][j / prime[i]];
		}
	}
	product[0] = 1;
	for (register int i = 1; i <= Q; i++){
		product[i] = product[i - 1] * prime[i];
	}
}
 
inline int get_id(ll n, ll m){
	return n <= sqrt_n ? id1[n] : id2[m / n];
}
 
inline ll cbrt(ll n){
	ll ans = cbrt((double)n);
	while (ans * ans * ans <= n) ans++;
	return ans - 1;
}
 
ll pi(ll n, ll m);
 
inline ll sum1(ll n){
	return n * (n + 1) / 2;
}
 
inline ll p2(ll n, ll m, ll k){
	ll a = pi(sqrt(n), k), ans = sum1(m - 1) - sum1(a - 1);
	for (register ll i = m + 1; i <= a; i++){
		ans += pi(n * inv[i], k);
	}
	return ans;
}
 
ll phi(ll n, ll m, ll k){
	if (n <= P && m <= K) return pre_phi[m][n];
	if (m <= Q) return n / product[m] * pre_phi[m][product[m]] + pre_phi[m][n % product[m]];
	if (n < prime[m] * prime[m]) return pi(n, k) - m + 1;
	if (n < prime[m] * prime[m] * prime[m]) return pi(n, k) + p2(n, m, k) - m + 1;
	return phi(n, m - 1, k) - phi(n * inv[m], m - 1, k);
}
 
inline ll p3(ll n, ll m, ll k){
	ll a = pi(cbrt(n), k), ans = 0;
	for (register ll i = m + 1; i <= a; i++){
		ans += p2(n * inv[i], i - 1, k);
	}
	return ans;
}
 
ll pi(ll n, ll m){
	if (n <= limit) return pi_list[n];
	int id = get_id(n, m);
	if (vis[id]) return big_pi_list[id];
	ll a = pi(sqrt(sqrt(n)), m);
	vis[id] = true;
	return big_pi_list[id] = phi(n, a, m) - p2(n, a, m) - p3(n, a, m) + a - 1;
}
inline int cal(int l,int r)
{
	if(l>r) return 0;
	return pi(r,100000000000)-pi(l-1,100000000000);
}
signed main(){
	ios::sync_with_stdio(false);
	cin.tie(0);
	init(1e11);
	int n;
	cin >> n;
	int ans=0;
	for(int i=1;prime[i]*prime[i]*prime[i]<=n;i++)
		++ans;
	for(int i=1;prime[i]*prime[i]<n;i++)
	{
	//	cout << prime[i] << "\n";
		ans+=cal(prime[i]+1,n/prime[i]);
	}
	cout << ans;
	return 0;
}
//https://codeforces.com/problemset/problem/665/F

