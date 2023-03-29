#include <iostream>
#include <bits/stdc++.h>
#include <cmath>
#include <vector>
#define ll long long int
#define mp make_pair
#define pb push_back
#define vi vector<int>
using namespace std;
const ll p =1e9 +7;
void solve(){
    ll n;
    cin>>n;
    ll ans=0;
    auto sumtill= [&](ll x){
       x%=p;
       ll sum=x*(x+1);
       sum%=p;
       sum*=(p+1)/2;
       sum%=p;
       return sum;
    };
    for(ll i=1;i*i<n;i++){
        if(i==n/i){
            continue;
        }
        ans+= n/i * i;
        ans%=p;
    }
    for(ll i=1;i*i<=n;i++){
        ans+=(sumtill(n/i) - sumtill(n/(i+1)))*i;
        ans%=p;
    }
    cout<<ans<<'\n';
}
int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    solve();
}
