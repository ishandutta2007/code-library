#include<bits/stdc++.h>
using namespace std;
typedef long long int lli;


lli exponentiation(lli a, lli b, lli c){
    if(a == 0)
        return 0;
    if(b == 0)
        return 1;
    lli ans ;
    if(b%2 == 0){
        lli smallAns = exponentiation(a, b/2, c);
        ans = (smallAns*smallAns) % c;
    }
    else{
        lli smallAns = exponentiation(a, b-1, c);
        ans = (a % c);
        ans = (ans*smallAns) % c;
    } 
    return (ans + c) % c;  //ans < 0 ? (ans + c)%c : ans;
}

lli wilson_Boring_factorials(lli n, lli p){
    
    if(n >= p)
        return 0; //becz n! will contain p in it and mod p will make it 0
    
    lli ans, fact = 1;
    //according to wilson thrm
    //n! modp = -[(n+1)^-1*.....(p-1)^-1] mod p = x modp , p>n
    //=>LHs = -x^-1 modp ------eqn1
    //applying fermet
    //LHS: -x^-1 modp = -x^(p-2) modp + p
    for(lli i=n+1; i<p; i++){
        fact = (fact * i) % p;
    }
    ans = exponentiation(fact, p-2, p); // we need fact^(p-2) mod p
    return p-ans;
}


int main(){
    
    int t;
    cin >> t;
    while(t--){
        lli n, p;
        cin >> n >> p;
        cout << wilson_Boring_factorials(n, p) << endl;
    }
    return 0;
}

//https://www.spoj.com/problems/BORING/

// NOT fast for https://www.spoj.com/problems/FACTMODP/